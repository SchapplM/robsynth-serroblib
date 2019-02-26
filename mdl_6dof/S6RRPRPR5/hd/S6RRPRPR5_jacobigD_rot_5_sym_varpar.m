% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRPR5_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:09
% EndTime: 2019-02-26 21:40:09
% DurationCPUTime: 0.04s
% Computational Cost: add. (15->8), mult. (54->22), div. (0->0), fcn. (58->8), ass. (0->16)
t170 = sin(pkin(11));
t172 = cos(pkin(11));
t174 = sin(qJ(2));
t176 = cos(qJ(2));
t179 = t176 * t170 + t174 * t172;
t181 = qJD(2) * t179;
t171 = sin(pkin(6));
t180 = qJD(1) * t171;
t178 = t170 * t174 - t172 * t176;
t177 = cos(qJ(1));
t175 = sin(qJ(1));
t173 = cos(pkin(6));
t168 = t178 * qJD(2);
t167 = t178 * t173;
t166 = t173 * t181;
t1 = [0, t177 * t180, 0, -t175 * t166 - t177 * t168 + (-t167 * t177 - t175 * t179) * qJD(1), 0, 0; 0, t175 * t180, 0, t177 * t166 - t175 * t168 + (-t167 * t175 + t177 * t179) * qJD(1), 0, 0; 0, 0, 0, t171 * t181, 0, 0;];
JgD_rot  = t1;
