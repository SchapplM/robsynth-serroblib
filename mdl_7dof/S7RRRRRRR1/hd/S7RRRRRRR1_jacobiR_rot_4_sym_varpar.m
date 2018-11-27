% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JR_rot [9x7]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S7RRRRRRR1_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiR_rot_4_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiR_rot_4_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:53
% EndTime: 2018-11-26 21:20:53
% DurationCPUTime: 0.06s
% Computational Cost: add. (37->22), mult. (118->39), div. (0->0), fcn. (180->8), ass. (0->30)
t90 = sin(qJ(4));
t92 = sin(qJ(2));
t112 = t92 * t90;
t94 = cos(qJ(4));
t111 = t92 * t94;
t95 = cos(qJ(3));
t110 = t92 * t95;
t97 = cos(qJ(1));
t109 = t92 * t97;
t91 = sin(qJ(3));
t93 = sin(qJ(1));
t108 = t93 * t91;
t107 = t93 * t95;
t96 = cos(qJ(2));
t106 = t96 * t90;
t105 = t96 * t91;
t104 = t96 * t94;
t103 = t97 * t91;
t102 = t97 * t95;
t86 = t96 * t107 + t103;
t101 = t93 * t111 - t86 * t90;
t100 = -t93 * t112 - t86 * t94;
t99 = -t94 * t110 + t106;
t98 = t90 * t110 + t104;
t88 = t96 * t102 - t108;
t87 = -t96 * t103 - t107;
t85 = t93 * t105 - t102;
t84 = t90 * t109 + t88 * t94;
t83 = t94 * t109 - t88 * t90;
t1 = [t100, t99 * t97, t87 * t94, t83, 0, 0, 0; t84, t99 * t93, -t85 * t94, t101, 0, 0, 0; 0, t95 * t104 + t112, -t91 * t111, -t98, 0, 0, 0; -t101, t98 * t97, -t87 * t90, -t84, 0, 0, 0; t83, t98 * t93, t85 * t90, t100, 0, 0, 0; 0, -t95 * t106 + t111, t91 * t112, t99, 0, 0, 0; t85, t92 * t103, -t88, 0, 0, 0, 0; t87, t92 * t108, -t86, 0, 0, 0, 0; 0, -t105, -t110, 0, 0, 0, 0;];
JR_rot  = t1;
