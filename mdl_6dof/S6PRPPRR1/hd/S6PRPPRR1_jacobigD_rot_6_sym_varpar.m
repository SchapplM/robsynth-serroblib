% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPPRR1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:35
% DurationCPUTime: 0.06s
% Computational Cost: add. (36->15), mult. (97->37), div. (0->0), fcn. (104->10), ass. (0->20)
t182 = sin(pkin(11));
t185 = cos(pkin(11));
t188 = sin(qJ(2));
t189 = cos(qJ(2));
t191 = t188 * t182 - t189 * t185;
t194 = t191 * qJD(2);
t192 = t182 * t189 + t185 * t188;
t177 = t192 * qJD(2);
t181 = pkin(12) + qJ(5);
t179 = sin(t181);
t184 = sin(pkin(6));
t193 = t184 * t179;
t187 = cos(pkin(6));
t186 = cos(pkin(10));
t183 = sin(pkin(10));
t180 = cos(t181);
t175 = t192 * t187;
t174 = t187 * t194;
t173 = t187 * t177;
t1 = [0, 0, 0, 0, -t183 * t173 - t186 * t194 (t183 * t174 - t186 * t177) * t179 + ((-t183 * t175 - t186 * t191) * t180 + t183 * t193) * qJD(5); 0, 0, 0, 0, t186 * t173 - t183 * t194 (-t186 * t174 - t183 * t177) * t179 + ((t186 * t175 - t183 * t191) * t180 - t186 * t193) * qJD(5); 0, 0, 0, 0, t184 * t177, t187 * qJD(5) * t179 + (t192 * qJD(5) * t180 - t179 * t194) * t184;];
JgD_rot  = t1;
