% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPPR5_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:51
% EndTime: 2019-02-26 22:05:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
t170 = sin(pkin(6));
t173 = sin(qJ(1));
t188 = t170 * t173;
t175 = cos(qJ(1));
t187 = t170 * t175;
t172 = sin(qJ(2));
t186 = t172 * t173;
t185 = t172 * t175;
t174 = cos(qJ(2));
t184 = t173 * t174;
t183 = t175 * t174;
t182 = qJD(1) * t170;
t169 = qJ(3) + pkin(11);
t167 = sin(t169);
t181 = qJD(2) * t167;
t180 = qJD(2) * t170;
t171 = cos(pkin(6));
t179 = t171 * t183 - t186;
t178 = t171 * t184 + t185;
t177 = t171 * t185 + t184;
t176 = -t171 * t186 + t183;
t168 = cos(t169);
t1 = [0, t175 * t182, t179 * qJD(1) + t176 * qJD(2), 0, 0 (t167 * t188 + t176 * t168) * qJD(3) - t178 * t181 + (-t177 * t167 - t168 * t187) * qJD(1); 0, t173 * t182, t178 * qJD(1) + t177 * qJD(2), 0, 0 (-t167 * t187 + t177 * t168) * qJD(3) + t179 * t181 + (t176 * t167 - t168 * t188) * qJD(1); 0, 0, t172 * t180, 0, 0, t174 * t167 * t180 + (t168 * t170 * t172 + t167 * t171) * qJD(3);];
JgD_rot  = t1;
