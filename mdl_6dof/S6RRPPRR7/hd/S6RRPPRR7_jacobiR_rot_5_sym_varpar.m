% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:32:07
% EndTime: 2019-02-26 21:32:07
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->18), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->25)
t82 = sin(pkin(6));
t84 = sin(qJ(5));
t99 = t82 * t84;
t87 = cos(qJ(5));
t98 = t82 * t87;
t88 = cos(qJ(2));
t97 = t82 * t88;
t89 = cos(qJ(1));
t96 = t82 * t89;
t85 = sin(qJ(2));
t86 = sin(qJ(1));
t95 = t86 * t85;
t94 = t86 * t88;
t93 = t89 * t85;
t92 = t89 * t88;
t83 = cos(pkin(6));
t77 = -t83 * t92 + t95;
t91 = -t77 * t84 + t87 * t96;
t90 = -t77 * t87 - t84 * t96;
t80 = -t83 * t95 + t92;
t79 = t83 * t94 + t93;
t78 = t83 * t93 + t94;
t76 = t79 * t87 - t86 * t99;
t75 = -t79 * t84 - t86 * t98;
t1 = [t90, t80 * t87, 0, 0, t75, 0; t76, t78 * t87, 0, 0, t91, 0; 0, t85 * t98, 0, 0, -t83 * t87 + t84 * t97, 0; -t91, -t80 * t84, 0, 0, -t76, 0; t75, -t78 * t84, 0, 0, t90, 0; 0, -t85 * t99, 0, 0, t83 * t84 + t87 * t97, 0; -t78, -t79, 0, 0, 0, 0; t80, -t77, 0, 0, 0, 0; 0, t97, 0, 0, 0, 0;];
JR_rot  = t1;
