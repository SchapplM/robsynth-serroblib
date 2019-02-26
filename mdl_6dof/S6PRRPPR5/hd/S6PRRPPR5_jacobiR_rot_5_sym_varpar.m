% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:42
% EndTime: 2019-02-26 20:00:42
% DurationCPUTime: 0.05s
% Computational Cost: add. (31->16), mult. (107->48), div. (0->0), fcn. (156->10), ass. (0->27)
t100 = sin(qJ(3));
t94 = sin(pkin(11));
t112 = t100 * t94;
t96 = sin(pkin(6));
t111 = t100 * t96;
t97 = cos(pkin(11));
t110 = t100 * t97;
t102 = cos(qJ(3));
t109 = t102 * t96;
t101 = sin(qJ(2));
t95 = sin(pkin(10));
t108 = t95 * t101;
t103 = cos(qJ(2));
t107 = t95 * t103;
t98 = cos(pkin(10));
t106 = t98 * t101;
t105 = t98 * t103;
t104 = t100 * t103;
t99 = cos(pkin(6));
t93 = t99 * t100 + t101 * t109;
t92 = -t99 * t108 + t105;
t91 = -t99 * t107 - t106;
t90 = t99 * t106 + t107;
t89 = t99 * t105 - t108;
t88 = t92 * t102 + t95 * t111;
t87 = t90 * t102 - t98 * t111;
t1 = [0, t91 * t112 + t92 * t97, t88 * t94, 0, 0, 0; 0, t89 * t112 + t90 * t97, t87 * t94, 0, 0, 0; 0 (t101 * t97 + t94 * t104) * t96, t93 * t94, 0, 0, 0; 0, t91 * t110 - t92 * t94, t88 * t97, 0, 0, 0; 0, t89 * t110 - t90 * t94, t87 * t97, 0, 0, 0; 0 (-t101 * t94 + t97 * t104) * t96, t93 * t97, 0, 0, 0; 0, t91 * t102, -t92 * t100 + t95 * t109, 0, 0, 0; 0, t89 * t102, -t90 * t100 - t98 * t109, 0, 0, 0; 0, t103 * t109, -t101 * t111 + t99 * t102, 0, 0, 0;];
JR_rot  = t1;
