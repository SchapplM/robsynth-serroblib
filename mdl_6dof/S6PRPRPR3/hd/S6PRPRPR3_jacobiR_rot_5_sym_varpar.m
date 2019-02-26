% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:44
% EndTime: 2019-02-26 19:47:44
% DurationCPUTime: 0.08s
% Computational Cost: add. (44->15), mult. (122->36), div. (0->0), fcn. (178->10), ass. (0->23)
t114 = sin(pkin(6));
t118 = sin(qJ(4));
t124 = t114 * t118;
t120 = cos(qJ(4));
t123 = t114 * t120;
t117 = cos(pkin(6));
t112 = sin(pkin(11));
t115 = cos(pkin(11));
t119 = sin(qJ(2));
t121 = cos(qJ(2));
t122 = t121 * t112 + t119 * t115;
t108 = t122 * t117;
t109 = t119 * t112 - t121 * t115;
t113 = sin(pkin(10));
t116 = cos(pkin(10));
t102 = t116 * t108 - t113 * t109;
t104 = -t113 * t108 - t116 * t109;
t107 = t109 * t117;
t106 = t122 * t114;
t105 = t109 * t114;
t103 = t113 * t107 - t116 * t122;
t101 = -t116 * t107 - t113 * t122;
t1 = [0, t104, 0, 0, 0, 0; 0, t102, 0, 0, 0, 0; 0, t106, 0, 0, 0, 0; 0, -t103 * t120, 0, t104 * t118 - t113 * t123, 0, 0; 0, -t101 * t120, 0, t102 * t118 + t116 * t123, 0, 0; 0, t105 * t120, 0, t106 * t118 - t117 * t120, 0, 0; 0, t103 * t118, 0, t104 * t120 + t113 * t124, 0, 0; 0, t101 * t118, 0, t102 * t120 - t116 * t124, 0, 0; 0, -t105 * t118, 0, t106 * t120 + t117 * t118, 0, 0;];
JR_rot  = t1;
