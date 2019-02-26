% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR13_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:01:09
% EndTime: 2019-02-26 22:01:09
% DurationCPUTime: 0.06s
% Computational Cost: add. (38->16), mult. (130->39), div. (0->0), fcn. (134->8), ass. (0->25)
t179 = sin(pkin(6));
t183 = sin(qJ(1));
t199 = t179 * t183;
t185 = cos(qJ(2));
t198 = t179 * t185;
t186 = cos(qJ(1));
t197 = t179 * t186;
t182 = sin(qJ(2));
t196 = t183 * t182;
t195 = t183 * t185;
t194 = t185 * t186;
t193 = t186 * t182;
t192 = qJD(1) * t179;
t184 = cos(qJ(4));
t191 = qJD(2) * t184;
t180 = cos(pkin(6));
t190 = t180 * t194 - t196;
t189 = t180 * t195 + t193;
t188 = t180 * t193 + t195;
t187 = -t180 * t196 + t194;
t181 = sin(qJ(4));
t178 = -t179 * t182 * t191 + (t180 * t184 - t181 * t198) * qJD(4);
t177 = (t189 * t181 + t184 * t199) * qJD(4) - t187 * t191 + (t181 * t197 - t190 * t184) * qJD(1);
t176 = (-t190 * t181 - t184 * t197) * qJD(4) - t188 * t191 + (t181 * t199 - t189 * t184) * qJD(1);
t1 = [0, t186 * t192, 0, -t188 * qJD(1) - t189 * qJD(2), t177, t177; 0, t183 * t192, 0, t187 * qJD(1) + t190 * qJD(2), t176, t176; 0, 0, 0, qJD(2) * t198, t178, t178;];
JgD_rot  = t1;
