% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP3_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:47:20
% EndTime: 2019-02-26 21:47:20
% DurationCPUTime: 0.03s
% Computational Cost: add. (35->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
t74 = sin(qJ(4));
t75 = sin(qJ(1));
t81 = t75 * t74;
t76 = cos(qJ(4));
t80 = t75 * t76;
t77 = cos(qJ(1));
t79 = t77 * t74;
t78 = t77 * t76;
t73 = qJ(2) + pkin(10);
t72 = cos(t73);
t71 = sin(t73);
t70 = t72 * t78 + t81;
t69 = -t72 * t79 + t80;
t68 = -t72 * t80 + t79;
t67 = t72 * t81 + t78;
t1 = [t68, -t71 * t78, 0, t69, 0, 0; t70, -t71 * t80, 0, -t67, 0, 0; 0, t72 * t76, 0, -t71 * t74, 0, 0; t67, t71 * t79, 0, -t70, 0, 0; t69, t71 * t81, 0, t68, 0, 0; 0, -t72 * t74, 0, -t71 * t76, 0, 0; -t75 * t71, t77 * t72, 0, 0, 0, 0; t77 * t71, t75 * t72, 0, 0, 0, 0; 0, t71, 0, 0, 0, 0;];
JR_rot  = t1;
