% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:43
% EndTime: 2019-02-26 20:11:44
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->12), mult. (60->30), div. (0->0), fcn. (62->8), ass. (0->21)
t177 = qJ(3) + qJ(4);
t175 = cos(t177);
t179 = sin(pkin(6));
t190 = t179 * t175;
t181 = cos(pkin(6));
t182 = sin(qJ(2));
t189 = t181 * t182;
t183 = cos(qJ(2));
t188 = t181 * t183;
t187 = qJD(2) * t175;
t186 = qJD(2) * t179;
t178 = sin(pkin(11));
t180 = cos(pkin(11));
t185 = t178 * t183 + t180 * t189;
t184 = -t178 * t189 + t180 * t183;
t176 = qJD(3) + qJD(4);
t174 = sin(t177);
t173 = t182 * t186;
t172 = t184 * qJD(2);
t171 = t185 * qJD(2);
t1 = [0, 0, t172, t172, 0 (-t184 * t174 + t178 * t190) * t176 + (-t178 * t188 - t180 * t182) * t187; 0, 0, t171, t171, 0 (-t185 * t174 - t180 * t190) * t176 + (-t178 * t182 + t180 * t188) * t187; 0, 0, t173, t173, 0, -t179 * t182 * t176 * t174 + (t176 * t181 + t183 * t186) * t175;];
JgD_rot  = t1;
