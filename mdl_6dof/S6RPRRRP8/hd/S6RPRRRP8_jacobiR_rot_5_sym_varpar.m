% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP8_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:11:54
% EndTime: 2019-02-26 21:11:54
% DurationCPUTime: 0.04s
% Computational Cost: add. (50->19), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t94 = qJ(3) + qJ(4);
t92 = sin(t94);
t97 = cos(qJ(5));
t106 = t92 * t97;
t95 = sin(qJ(5));
t96 = sin(qJ(1));
t105 = t96 * t95;
t104 = t96 * t97;
t98 = cos(qJ(1));
t103 = t98 * t92;
t102 = t98 * t95;
t101 = t98 * t97;
t93 = cos(t94);
t100 = t93 * t105;
t99 = t93 * t101;
t91 = t96 * t92;
t90 = t92 * t95;
t89 = t93 * t102;
t88 = t93 * t104;
t87 = t92 * t101 - t105;
t86 = t92 * t102 + t104;
t85 = t92 * t104 + t102;
t84 = -t92 * t105 + t101;
t1 = [t87, 0, t88, t88, t84, 0; t85, 0, -t99, -t99, t86, 0; 0, 0, -t106, -t106, -t93 * t95, 0; -t86, 0, -t100, -t100, -t85, 0; t84, 0, t89, t89, t87, 0; 0, 0, t90, t90, -t93 * t97, 0; -t98 * t93, 0, t91, t91, 0, 0; -t96 * t93, 0, -t103, -t103, 0, 0; 0, 0, t93, t93, 0, 0;];
JR_rot  = t1;
