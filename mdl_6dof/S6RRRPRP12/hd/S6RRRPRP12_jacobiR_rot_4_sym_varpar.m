% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP12_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:15:27
% EndTime: 2019-02-26 22:15:27
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->15), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->25)
t101 = sin(qJ(1));
t97 = sin(pkin(6));
t114 = t101 * t97;
t103 = cos(qJ(2));
t113 = t103 * t97;
t104 = cos(qJ(1));
t112 = t104 * t97;
t100 = sin(qJ(2));
t111 = t97 * t100;
t110 = t101 * t100;
t109 = t101 * t103;
t108 = t104 * t100;
t107 = t104 * t103;
t102 = cos(qJ(3));
t98 = cos(pkin(6));
t94 = t98 * t108 + t109;
t99 = sin(qJ(3));
t106 = t94 * t102 - t99 * t112;
t105 = t102 * t112 + t94 * t99;
t96 = -t98 * t110 + t107;
t95 = t98 * t109 + t108;
t93 = t98 * t107 - t110;
t92 = t96 * t102 + t99 * t114;
t91 = -t102 * t114 + t96 * t99;
t1 = [t93, t96, 0, 0, 0, 0; t95, t94, 0, 0, 0, 0; 0, t111, 0, 0, 0, 0; t106, t95 * t102, t91, 0, 0, 0; -t92, -t93 * t102, t105, 0, 0, 0; 0, -t102 * t113, -t98 * t102 + t99 * t111, 0, 0, 0; -t105, -t95 * t99, t92, 0, 0, 0; t91, t93 * t99, t106, 0, 0, 0; 0, t99 * t113, t102 * t111 + t98 * t99, 0, 0, 0;];
JR_rot  = t1;
