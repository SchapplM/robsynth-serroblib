% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRPR3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:39
% EndTime: 2019-02-26 19:47:39
% DurationCPUTime: 0.06s
% Computational Cost: add. (27->14), mult. (97->37), div. (0->0), fcn. (104->10), ass. (0->19)
t175 = sin(pkin(11));
t178 = cos(pkin(11));
t182 = sin(qJ(2));
t184 = cos(qJ(2));
t186 = t182 * t175 - t184 * t178;
t189 = t186 * qJD(2);
t187 = t175 * t184 + t178 * t182;
t173 = t187 * qJD(2);
t177 = sin(pkin(6));
t183 = cos(qJ(4));
t188 = t177 * t183;
t181 = sin(qJ(4));
t180 = cos(pkin(6));
t179 = cos(pkin(10));
t176 = sin(pkin(10));
t171 = t187 * t180;
t170 = t180 * t189;
t169 = t180 * t173;
t1 = [0, 0, 0, -t176 * t169 - t179 * t189, 0 (t176 * t170 - t179 * t173) * t183 + (-(-t176 * t171 - t179 * t186) * t181 + t176 * t188) * qJD(4); 0, 0, 0, t179 * t169 - t176 * t189, 0 (-t179 * t170 - t176 * t173) * t183 + (-(t179 * t171 - t176 * t186) * t181 - t179 * t188) * qJD(4); 0, 0, 0, t177 * t173, 0, t180 * qJD(4) * t183 + (-t187 * qJD(4) * t181 - t183 * t189) * t177;];
JgD_rot  = t1;
