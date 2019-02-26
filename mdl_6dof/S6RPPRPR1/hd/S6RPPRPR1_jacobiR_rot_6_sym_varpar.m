% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRPR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:25:29
% EndTime: 2019-02-26 20:25:30
% DurationCPUTime: 0.04s
% Computational Cost: add. (83->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->21)
t80 = pkin(10) + qJ(4);
t74 = sin(t80);
t81 = qJ(1) + pkin(9);
t75 = sin(t81);
t88 = t75 * t74;
t79 = pkin(11) + qJ(6);
t76 = cos(t79);
t87 = t75 * t76;
t73 = sin(t79);
t77 = cos(t80);
t86 = t77 * t73;
t85 = t77 * t76;
t78 = cos(t81);
t84 = t78 * t74;
t83 = t78 * t76;
t82 = t78 * t77;
t72 = t75 * t73 + t76 * t82;
t71 = -t73 * t82 + t87;
t70 = t78 * t73 - t75 * t85;
t69 = t75 * t86 + t83;
t1 = [t70, 0, 0, -t74 * t83, 0, t71; t72, 0, 0, -t74 * t87, 0, -t69; 0, 0, 0, t85, 0, -t74 * t73; t69, 0, 0, t73 * t84, 0, -t72; t71, 0, 0, t73 * t88, 0, t70; 0, 0, 0, -t86, 0, -t74 * t76; -t88, 0, 0, t82, 0, 0; t84, 0, 0, t75 * t77, 0, 0; 0, 0, 0, t74, 0, 0;];
JR_rot  = t1;
