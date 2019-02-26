% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:48:47
% EndTime: 2019-02-26 22:48:47
% DurationCPUTime: 0.09s
% Computational Cost: add. (194->25), mult. (82->20), div. (0->0), fcn. (141->6), ass. (0->19)
t92 = qJ(3) + qJ(4) + qJ(5) + qJ(6);
t90 = sin(t92);
t93 = sin(qJ(2));
t103 = t93 * t90;
t91 = cos(t92);
t102 = t93 * t91;
t94 = sin(qJ(1));
t101 = t94 * t93;
t95 = cos(qJ(2));
t100 = t95 * t90;
t99 = t95 * t91;
t96 = cos(qJ(1));
t98 = t96 * t93;
t97 = t96 * t95;
t89 = t94 * t90 + t91 * t97;
t88 = -t90 * t97 + t94 * t91;
t87 = t96 * t90 - t94 * t99;
t86 = t100 * t94 + t96 * t91;
t1 = [t87, -t91 * t98, t88, t88, t88, t88; t89, -t91 * t101, -t86, -t86, -t86, -t86; 0, t99, -t103, -t103, -t103, -t103; t86, t90 * t98, -t89, -t89, -t89, -t89; t88, t90 * t101, t87, t87, t87, t87; 0, -t100, -t102, -t102, -t102, -t102; -t101, t97, 0, 0, 0, 0; t98, t94 * t95, 0, 0, 0, 0; 0, t93, 0, 0, 0, 0;];
JR_rot  = t1;
