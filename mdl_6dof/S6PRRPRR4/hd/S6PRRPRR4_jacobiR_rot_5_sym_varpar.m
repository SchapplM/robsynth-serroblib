% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:52
% EndTime: 2019-02-26 20:05:52
% DurationCPUTime: 0.06s
% Computational Cost: add. (69->27), mult. (203->47), div. (0->0), fcn. (292->10), ass. (0->36)
t94 = sin(pkin(6));
t98 = sin(qJ(3));
t117 = t94 * t98;
t99 = sin(qJ(2));
t116 = t94 * t99;
t96 = cos(pkin(6));
t115 = t96 * t99;
t101 = cos(qJ(3));
t114 = t101 * t94;
t102 = cos(qJ(2));
t113 = t102 * t94;
t93 = sin(pkin(11));
t112 = t93 * t102;
t95 = cos(pkin(11));
t111 = t95 * t102;
t100 = cos(qJ(5));
t86 = t95 * t115 + t112;
t81 = t95 * t114 + t86 * t98;
t82 = t86 * t101 - t95 * t117;
t97 = sin(qJ(5));
t110 = t82 * t100 + t81 * t97;
t109 = t81 * t100 - t82 * t97;
t88 = -t93 * t115 + t111;
t83 = -t93 * t114 + t88 * t98;
t84 = t88 * t101 + t93 * t117;
t108 = t84 * t100 + t83 * t97;
t107 = t83 * t100 - t84 * t97;
t89 = -t96 * t101 + t98 * t116;
t90 = t99 * t114 + t96 * t98;
t106 = t90 * t100 + t89 * t97;
t105 = t89 * t100 - t90 * t97;
t104 = t100 * t98 - t101 * t97;
t103 = t100 * t101 + t97 * t98;
t87 = -t96 * t112 - t95 * t99;
t85 = t96 * t111 - t93 * t99;
t1 = [0, t103 * t87, -t107, 0, t107, 0; 0, t103 * t85, -t109, 0, t109, 0; 0, t103 * t113, -t105, 0, t105, 0; 0, t104 * t87, t108, 0, -t108, 0; 0, t104 * t85, t110, 0, -t110, 0; 0, t104 * t113, t106, 0, -t106, 0; 0, -t88, 0, 0, 0, 0; 0, -t86, 0, 0, 0, 0; 0, -t116, 0, 0, 0, 0;];
JR_rot  = t1;
