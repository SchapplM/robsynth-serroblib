% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JgD_rot = S6PRPRRP1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:50:18
% EndTime: 2019-02-26 19:50:19
% DurationCPUTime: 0.06s
% Computational Cost: add. (27->14), mult. (97->37), div. (0->0), fcn. (104->10), ass. (0->19)
t183 = sin(pkin(11));
t186 = cos(pkin(11));
t190 = sin(qJ(2));
t192 = cos(qJ(2));
t194 = t190 * t183 - t192 * t186;
t197 = t194 * qJD(2);
t195 = t183 * t192 + t186 * t190;
t181 = t195 * qJD(2);
t185 = sin(pkin(6));
t189 = sin(qJ(4));
t196 = t185 * t189;
t191 = cos(qJ(4));
t188 = cos(pkin(6));
t187 = cos(pkin(10));
t184 = sin(pkin(10));
t179 = t195 * t188;
t178 = t188 * t197;
t177 = t188 * t181;
t1 = [0, 0, 0, -t184 * t177 - t187 * t197 (t184 * t178 - t187 * t181) * t189 + ((-t184 * t179 - t187 * t194) * t191 + t184 * t196) * qJD(4), 0; 0, 0, 0, t187 * t177 - t184 * t197 (-t187 * t178 - t184 * t181) * t189 + ((t187 * t179 - t184 * t194) * t191 - t187 * t196) * qJD(4), 0; 0, 0, 0, t185 * t181, t188 * qJD(4) * t189 + (t195 * qJD(4) * t191 - t189 * t197) * t185, 0;];
JgD_rot  = t1;
