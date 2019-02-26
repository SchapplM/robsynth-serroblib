% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR8_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:39
% EndTime: 2019-02-26 22:07:39
% DurationCPUTime: 0.04s
% Computational Cost: add. (30->18), mult. (87->31), div. (0->0), fcn. (134->8), ass. (0->25)
t94 = sin(pkin(6));
t97 = sin(qJ(2));
t111 = t94 * t97;
t99 = cos(qJ(3));
t110 = t94 * t99;
t98 = sin(qJ(1));
t109 = t98 * t97;
t100 = cos(qJ(2));
t108 = t100 * t94;
t101 = cos(qJ(1));
t107 = t101 * t94;
t106 = t101 * t97;
t105 = t98 * t100;
t104 = t101 * t100;
t95 = cos(pkin(6));
t90 = t95 * t106 + t105;
t96 = sin(qJ(3));
t103 = t99 * t107 + t90 * t96;
t102 = -t96 * t107 + t90 * t99;
t92 = -t95 * t109 + t104;
t91 = -t95 * t105 - t106;
t89 = -t95 * t104 + t109;
t88 = t98 * t94 * t96 + t92 * t99;
t87 = -t98 * t110 + t92 * t96;
t1 = [-t103, t91 * t96, t88, 0, 0, 0; t87, -t89 * t96, t102, 0, 0, 0; 0, t96 * t108, t97 * t110 + t95 * t96, 0, 0, 0; t102, -t91 * t99, t87, 0, 0, 0; -t88, t89 * t99, t103, 0, 0, 0; 0, -t99 * t108, t96 * t111 - t95 * t99, 0, 0, 0; t89, -t92, 0, 0, 0, 0; t91, -t90, 0, 0, 0, 0; 0, -t111, 0, 0, 0, 0;];
JR_rot  = t1;
