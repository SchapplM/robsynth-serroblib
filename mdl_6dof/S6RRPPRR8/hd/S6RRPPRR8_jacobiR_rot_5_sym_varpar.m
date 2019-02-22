% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:17
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR8_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:17:02
% EndTime: 2019-02-22 11:17:02
% DurationCPUTime: 0.06s
% Computational Cost: add. (36->21), mult. (108->32), div. (0->0), fcn. (161->8), ass. (0->23)
t79 = sin(qJ(2));
t76 = sin(pkin(10));
t77 = cos(pkin(10));
t78 = sin(qJ(5));
t81 = cos(qJ(5));
t87 = t76 * t81 - t77 * t78;
t85 = t87 * t79;
t80 = sin(qJ(1));
t82 = cos(qJ(2));
t91 = t80 * t82;
t83 = cos(qJ(1));
t90 = t83 * t82;
t71 = t76 * t91 + t83 * t77;
t72 = -t83 * t76 + t77 * t91;
t89 = t71 * t81 - t72 * t78;
t88 = -t71 * t78 - t72 * t81;
t86 = t76 * t78 + t77 * t81;
t84 = t86 * t79;
t74 = t80 * t76 + t77 * t90;
t73 = t76 * t90 - t80 * t77;
t70 = t73 * t78 + t74 * t81;
t69 = t73 * t81 - t74 * t78;
t1 = [t88, -t83 * t84, 0, 0, t69, 0; t70, -t80 * t84, 0, 0, t89, 0; 0, t86 * t82, 0, 0, t85, 0; -t89, -t83 * t85, 0, 0, -t70, 0; t69, -t80 * t85, 0, 0, t88, 0; 0, t87 * t82, 0, 0, -t84, 0; t80 * t79, -t90, 0, 0, 0, 0; -t83 * t79, -t91, 0, 0, 0, 0; 0, -t79, 0, 0, 0, 0;];
JR_rot  = t1;
