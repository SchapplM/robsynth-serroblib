% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:39
% EndTime: 2019-02-26 20:10:39
% DurationCPUTime: 0.04s
% Computational Cost: add. (38->12), mult. (60->30), div. (0->0), fcn. (62->8), ass. (0->21)
t174 = qJ(3) + qJ(4) + pkin(12);
t172 = sin(t174);
t177 = sin(pkin(6));
t188 = t177 * t172;
t179 = cos(pkin(6));
t180 = sin(qJ(2));
t187 = t179 * t180;
t181 = cos(qJ(2));
t186 = t179 * t181;
t185 = qJD(2) * t172;
t184 = qJD(2) * t177;
t176 = sin(pkin(11));
t178 = cos(pkin(11));
t183 = t176 * t181 + t178 * t187;
t182 = -t176 * t187 + t178 * t181;
t175 = qJD(3) + qJD(4);
t173 = cos(t174);
t171 = t180 * t184;
t170 = t182 * qJD(2);
t169 = t183 * qJD(2);
t1 = [0, 0, t170, t170, 0 (t182 * t173 + t176 * t188) * t175 + (-t176 * t186 - t178 * t180) * t185; 0, 0, t169, t169, 0 (t183 * t173 - t178 * t188) * t175 + (-t176 * t180 + t178 * t186) * t185; 0, 0, t171, t171, 0, t177 * t180 * t175 * t173 + (t175 * t179 + t181 * t184) * t172;];
JgD_rot  = t1;
