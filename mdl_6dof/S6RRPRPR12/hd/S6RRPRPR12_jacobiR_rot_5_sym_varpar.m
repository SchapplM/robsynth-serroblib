% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR12_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:12
% EndTime: 2019-02-26 21:44:12
% DurationCPUTime: 0.05s
% Computational Cost: add. (52->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
t95 = sin(pkin(6));
t97 = sin(qJ(2));
t110 = t95 * t97;
t98 = sin(qJ(1));
t109 = t95 * t98;
t99 = cos(qJ(2));
t108 = t95 * t99;
t107 = t98 * t97;
t106 = t98 * t99;
t100 = cos(qJ(1));
t105 = t100 * t95;
t104 = t100 * t97;
t103 = t100 * t99;
t96 = cos(pkin(6));
t86 = -t96 * t103 + t107;
t94 = qJ(4) + pkin(11);
t92 = sin(t94);
t93 = cos(t94);
t102 = t93 * t105 - t86 * t92;
t101 = t92 * t105 + t86 * t93;
t89 = -t96 * t107 + t103;
t88 = t96 * t106 + t104;
t87 = t96 * t104 + t106;
t85 = t93 * t109 + t88 * t92;
t84 = -t92 * t109 + t88 * t93;
t1 = [t102, t89 * t92, 0, t84, 0, 0; t85, t87 * t92, 0, t101, 0, 0; 0, t92 * t110, 0, -t93 * t108 - t96 * t92, 0, 0; -t101, t89 * t93, 0, -t85, 0, 0; t84, t87 * t93, 0, t102, 0, 0; 0, t93 * t110, 0, t92 * t108 - t96 * t93, 0, 0; -t87, -t88, 0, 0, 0, 0; t89, -t86, 0, 0, 0, 0; 0, t108, 0, 0, 0, 0;];
JR_rot  = t1;
