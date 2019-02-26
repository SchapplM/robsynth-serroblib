% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRPR2_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:04
% EndTime: 2019-02-26 19:47:04
% DurationCPUTime: 0.06s
% Computational Cost: add. (27->14), mult. (97->37), div. (0->0), fcn. (104->10), ass. (0->19)
t181 = sin(pkin(11));
t184 = cos(pkin(11));
t188 = sin(qJ(2));
t190 = cos(qJ(2));
t192 = t188 * t181 - t190 * t184;
t195 = t192 * qJD(2);
t193 = t181 * t190 + t184 * t188;
t179 = t193 * qJD(2);
t183 = sin(pkin(6));
t187 = sin(qJ(4));
t194 = t183 * t187;
t189 = cos(qJ(4));
t186 = cos(pkin(6));
t185 = cos(pkin(10));
t182 = sin(pkin(10));
t177 = t193 * t186;
t176 = t186 * t195;
t175 = t186 * t179;
t1 = [0, 0, 0, -t182 * t175 - t185 * t195, 0 (t182 * t176 - t185 * t179) * t187 + ((-t182 * t177 - t185 * t192) * t189 + t182 * t194) * qJD(4); 0, 0, 0, t185 * t175 - t182 * t195, 0 (-t185 * t176 - t182 * t179) * t187 + ((t185 * t177 - t182 * t192) * t189 - t185 * t194) * qJD(4); 0, 0, 0, t183 * t179, 0, t186 * qJD(4) * t187 + (t193 * qJD(4) * t189 - t187 * t195) * t183;];
JgD_rot  = t1;
