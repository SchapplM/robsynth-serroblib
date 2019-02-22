% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:22
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPPR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:22:12
% EndTime: 2019-02-22 10:22:12
% DurationCPUTime: 0.04s
% Computational Cost: add. (83->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->21)
t83 = qJ(3) + pkin(10);
t77 = sin(t83);
t84 = qJ(1) + pkin(9);
t78 = sin(t84);
t91 = t78 * t77;
t82 = pkin(11) + qJ(6);
t79 = cos(t82);
t90 = t78 * t79;
t76 = sin(t82);
t80 = cos(t83);
t89 = t80 * t76;
t88 = t80 * t79;
t81 = cos(t84);
t87 = t81 * t77;
t86 = t81 * t79;
t85 = t81 * t80;
t75 = t78 * t76 + t79 * t85;
t74 = -t76 * t85 + t90;
t73 = t81 * t76 - t78 * t88;
t72 = t78 * t89 + t86;
t1 = [t73, 0, -t77 * t86, 0, 0, t74; t75, 0, -t77 * t90, 0, 0, -t72; 0, 0, t88, 0, 0, -t77 * t76; t72, 0, t76 * t87, 0, 0, -t75; t74, 0, t76 * t91, 0, 0, t73; 0, 0, -t89, 0, 0, -t77 * t79; -t91, 0, t85, 0, 0, 0; t87, 0, t78 * t80, 0, 0, 0; 0, 0, t77, 0, 0, 0;];
JR_rot  = t1;
