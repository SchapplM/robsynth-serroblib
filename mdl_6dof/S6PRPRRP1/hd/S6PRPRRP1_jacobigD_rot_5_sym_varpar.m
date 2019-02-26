% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRP1_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_jacobigD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:50:18
% EndTime: 2019-02-26 19:50:19
% DurationCPUTime: 0.06s
% Computational Cost: add. (27->14), mult. (97->37), div. (0->0), fcn. (104->10), ass. (0->19)
t172 = sin(pkin(11));
t175 = cos(pkin(11));
t179 = sin(qJ(2));
t181 = cos(qJ(2));
t183 = t179 * t172 - t181 * t175;
t186 = t183 * qJD(2);
t184 = t172 * t181 + t175 * t179;
t170 = t184 * qJD(2);
t174 = sin(pkin(6));
t178 = sin(qJ(4));
t185 = t174 * t178;
t180 = cos(qJ(4));
t177 = cos(pkin(6));
t176 = cos(pkin(10));
t173 = sin(pkin(10));
t168 = t184 * t177;
t167 = t177 * t186;
t166 = t177 * t170;
t1 = [0, 0, 0, -t173 * t166 - t176 * t186 (t173 * t167 - t176 * t170) * t178 + ((-t173 * t168 - t176 * t183) * t180 + t173 * t185) * qJD(4), 0; 0, 0, 0, t176 * t166 - t173 * t186 (-t176 * t167 - t173 * t170) * t178 + ((t176 * t168 - t173 * t183) * t180 - t176 * t185) * qJD(4), 0; 0, 0, 0, t174 * t170, t177 * qJD(4) * t178 + (t184 * qJD(4) * t180 - t178 * t186) * t174, 0;];
JgD_rot  = t1;
