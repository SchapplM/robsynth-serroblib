% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRP1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:13
% EndTime: 2019-02-26 20:15:13
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->12), mult. (60->30), div. (0->0), fcn. (62->8), ass. (0->21)
t180 = qJ(3) + qJ(4);
t177 = sin(t180);
t182 = sin(pkin(6));
t193 = t182 * t177;
t184 = cos(pkin(6));
t185 = sin(qJ(2));
t192 = t184 * t185;
t186 = cos(qJ(2));
t191 = t184 * t186;
t190 = qJD(2) * t177;
t189 = qJD(2) * t182;
t181 = sin(pkin(11));
t183 = cos(pkin(11));
t188 = t181 * t186 + t183 * t192;
t187 = -t181 * t192 + t183 * t186;
t179 = qJD(3) + qJD(4);
t178 = cos(t180);
t176 = t185 * t189;
t175 = t187 * qJD(2);
t174 = t188 * qJD(2);
t1 = [0, 0, t175, t175 (t187 * t178 + t181 * t193) * t179 + (-t181 * t191 - t183 * t185) * t190, 0; 0, 0, t174, t174 (t188 * t178 - t183 * t193) * t179 + (-t181 * t185 + t183 * t191) * t190, 0; 0, 0, t176, t176, t182 * t185 * t179 * t178 + (t179 * t184 + t186 * t189) * t177, 0;];
JgD_rot  = t1;
