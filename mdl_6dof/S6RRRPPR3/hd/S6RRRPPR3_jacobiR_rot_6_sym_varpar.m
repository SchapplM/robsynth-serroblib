% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:04:32
% EndTime: 2019-02-26 22:04:32
% DurationCPUTime: 0.04s
% Computational Cost: add. (49->18), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t102 = qJ(2) + qJ(3);
t100 = sin(t102);
t103 = sin(qJ(6));
t115 = t100 * t103;
t104 = sin(qJ(1));
t114 = t104 * t100;
t113 = t104 * t103;
t105 = cos(qJ(6));
t112 = t104 * t105;
t106 = cos(qJ(1));
t111 = t106 * t100;
t110 = t106 * t103;
t109 = t106 * t105;
t101 = cos(t102);
t108 = t101 * t113;
t107 = t101 * t110;
t99 = t100 * t105;
t98 = t101 * t109;
t97 = t101 * t112;
t96 = t100 * t109 - t113;
t95 = -t100 * t110 - t112;
t94 = -t100 * t112 - t110;
t93 = t100 * t113 - t109;
t1 = [t94, t98, t98, 0, 0, t95; t96, t97, t97, 0, 0, -t93; 0, t99, t99, 0, 0, t101 * t103; t93, -t107, -t107, 0, 0, -t96; t95, -t108, -t108, 0, 0, t94; 0, -t115, -t115, 0, 0, t101 * t105; -t104 * t101, -t111, -t111, 0, 0, 0; t106 * t101, -t114, -t114, 0, 0, 0; 0, t101, t101, 0, 0, 0;];
JR_rot  = t1;
