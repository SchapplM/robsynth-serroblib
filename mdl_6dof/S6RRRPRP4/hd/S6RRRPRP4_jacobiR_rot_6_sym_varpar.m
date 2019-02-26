% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:05
% EndTime: 2019-02-26 22:11:05
% DurationCPUTime: 0.07s
% Computational Cost: add. (51->20), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t111 = qJ(2) + qJ(3);
t109 = sin(t111);
t114 = cos(qJ(5));
t124 = t109 * t114;
t113 = sin(qJ(1));
t123 = t113 * t109;
t112 = sin(qJ(5));
t122 = t113 * t112;
t121 = t113 * t114;
t115 = cos(qJ(1));
t120 = t115 * t109;
t119 = t115 * t112;
t118 = t115 * t114;
t110 = cos(t111);
t117 = t110 * t121;
t116 = t110 * t118;
t108 = t109 * t112;
t107 = t110 * t119;
t106 = t110 * t122;
t105 = -t109 * t122 + t118;
t104 = t109 * t121 + t119;
t103 = t109 * t119 + t121;
t102 = -t109 * t118 + t122;
t1 = [t105, t107, t107, 0, -t102, 0; t103, t106, t106, 0, t104, 0; 0, t108, t108, 0, -t110 * t114, 0; -t113 * t110, -t120, -t120, 0, 0, 0; t115 * t110, -t123, -t123, 0, 0, 0; 0, t110, t110, 0, 0, 0; t104, -t116, -t116, 0, t103, 0; t102, -t117, -t117, 0, -t105, 0; 0, -t124, -t124, 0, -t110 * t112, 0;];
JR_rot  = t1;
