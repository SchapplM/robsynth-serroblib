% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:53
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:53:29
% EndTime: 2019-02-22 11:53:29
% DurationCPUTime: 0.09s
% Computational Cost: add. (36->19), mult. (108->32), div. (0->0), fcn. (161->8), ass. (0->23)
t85 = sin(qJ(2));
t82 = sin(pkin(10));
t83 = cos(pkin(10));
t84 = sin(qJ(3));
t87 = cos(qJ(3));
t92 = t82 * t84 + t83 * t87;
t91 = t92 * t85;
t86 = sin(qJ(1));
t88 = cos(qJ(2));
t99 = t86 * t88;
t89 = cos(qJ(1));
t98 = t89 * t88;
t76 = -t84 * t99 - t89 * t87;
t77 = -t89 * t84 + t87 * t99;
t97 = t76 * t83 + t77 * t82;
t78 = t84 * t98 - t86 * t87;
t79 = t86 * t84 + t87 * t98;
t96 = t78 * t82 + t79 * t83;
t95 = t76 * t82 - t77 * t83;
t94 = t78 * t83 - t79 * t82;
t93 = t82 * t87 - t83 * t84;
t90 = t93 * t85;
t1 = [t95, -t89 * t91, -t94, 0, 0, 0; t96, -t86 * t91, t97, 0, 0, 0; 0, t92 * t88, t90, 0, 0, 0; t97, t89 * t90, t96, 0, 0, 0; t94, t86 * t90, -t95, 0, 0, 0; 0, -t93 * t88, t91, 0, 0, 0; t86 * t85, -t98, 0, 0, 0, 0; -t89 * t85, -t99, 0, 0, 0, 0; 0, -t85, 0, 0, 0, 0;];
JR_rot  = t1;
