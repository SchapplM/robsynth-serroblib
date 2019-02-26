% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPPRRR1_jacobig_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobig_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobig_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:53
% EndTime: 2019-02-26 19:38:54
% DurationCPUTime: 0.04s
% Computational Cost: add. (18->16), mult. (52->38), div. (0->0), fcn. (73->12), ass. (0->19)
t67 = sin(pkin(12));
t76 = cos(pkin(6));
t80 = t67 * t76;
t69 = sin(pkin(7));
t70 = sin(pkin(6));
t79 = t69 * t70;
t75 = cos(pkin(7));
t78 = t70 * t75;
t73 = cos(pkin(12));
t77 = t73 * t76;
t74 = cos(pkin(8));
t72 = cos(pkin(13));
t71 = cos(pkin(14));
t68 = sin(pkin(8));
t66 = sin(pkin(13));
t65 = sin(pkin(14));
t64 = -t73 * t66 - t72 * t80;
t63 = -t67 * t66 + t72 * t77;
t1 = [0, 0, 0 -(-(-t66 * t80 + t73 * t72) * t65 + (t64 * t75 + t67 * t79) * t71) * t68 + (-t64 * t69 + t67 * t78) * t74, 0, 0; 0, 0, 0 -(-(t66 * t77 + t67 * t72) * t65 + (t63 * t75 - t73 * t79) * t71) * t68 + (-t63 * t69 - t73 * t78) * t74, 0, 0; 0, 0, 0 -(t76 * t69 * t71 + (t71 * t72 * t75 - t65 * t66) * t70) * t68 + (-t72 * t79 + t76 * t75) * t74, 0, 0;];
Jg_rot  = t1;
