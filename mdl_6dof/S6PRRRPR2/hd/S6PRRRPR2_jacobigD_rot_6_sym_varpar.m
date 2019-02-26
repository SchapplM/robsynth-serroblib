% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR2
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
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR2_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:07
% EndTime: 2019-02-26 20:11:07
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->12), mult. (60->30), div. (0->0), fcn. (62->8), ass. (0->21)
t179 = qJ(3) + qJ(4);
t176 = sin(t179);
t181 = sin(pkin(6));
t192 = t181 * t176;
t183 = cos(pkin(6));
t184 = sin(qJ(2));
t191 = t183 * t184;
t185 = cos(qJ(2));
t190 = t183 * t185;
t189 = qJD(2) * t176;
t188 = qJD(2) * t181;
t180 = sin(pkin(11));
t182 = cos(pkin(11));
t187 = t180 * t185 + t182 * t191;
t186 = -t180 * t191 + t182 * t185;
t178 = qJD(3) + qJD(4);
t177 = cos(t179);
t175 = t184 * t188;
t174 = t186 * qJD(2);
t173 = t187 * qJD(2);
t1 = [0, 0, t174, t174, 0 (t186 * t177 + t180 * t192) * t178 + (-t180 * t190 - t182 * t184) * t189; 0, 0, t173, t173, 0 (t187 * t177 - t182 * t192) * t178 + (-t180 * t184 + t182 * t190) * t189; 0, 0, t175, t175, 0, t181 * t184 * t178 * t177 + (t178 * t183 + t185 * t188) * t176;];
JgD_rot  = t1;
