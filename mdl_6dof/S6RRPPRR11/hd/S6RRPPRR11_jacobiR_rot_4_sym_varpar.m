% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:19
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR11_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobiR_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:18:58
% EndTime: 2019-02-22 11:18:58
% DurationCPUTime: 0.04s
% Computational Cost: add. (16->10), mult. (57->26), div. (0->0), fcn. (88->8), ass. (0->20)
t67 = sin(pkin(6));
t70 = sin(qJ(2));
t80 = t67 * t70;
t71 = sin(qJ(1));
t79 = t67 * t71;
t73 = cos(qJ(1));
t78 = t67 * t73;
t77 = t71 * t70;
t72 = cos(qJ(2));
t76 = t71 * t72;
t75 = t73 * t70;
t74 = t73 * t72;
t69 = cos(pkin(6));
t68 = cos(pkin(11));
t66 = sin(pkin(11));
t65 = -t69 * t77 + t74;
t64 = t69 * t76 + t75;
t63 = t69 * t75 + t76;
t62 = t69 * t74 - t77;
t1 = [t62 * t66 + t68 * t78, t65 * t66, 0, 0, 0, 0; t64 * t66 + t68 * t79, t63 * t66, 0, 0, 0, 0; 0, t66 * t80, 0, 0, 0, 0; t62 * t68 - t66 * t78, t65 * t68, 0, 0, 0, 0; t64 * t68 - t66 * t79, t63 * t68, 0, 0, 0, 0; 0, t68 * t80, 0, 0, 0, 0; -t63, -t64, 0, 0, 0, 0; t65, t62, 0, 0, 0, 0; 0, t67 * t72, 0, 0, 0, 0;];
JR_rot  = t1;
