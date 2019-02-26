% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR9_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobia_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:23
% EndTime: 2019-02-26 20:53:23
% DurationCPUTime: 0.10s
% Computational Cost: add. (218->26), mult. (618->57), div. (26->9), fcn. (869->15), ass. (0->45)
t74 = sin(pkin(7));
t72 = sin(pkin(13));
t76 = cos(pkin(13));
t80 = sin(qJ(3));
t82 = cos(qJ(3));
t84 = t82 * t72 + t80 * t76;
t59 = t84 * t74;
t78 = cos(pkin(7));
t61 = t84 * t78;
t79 = cos(pkin(6));
t73 = sin(pkin(12));
t83 = cos(qJ(1));
t86 = t83 * t73;
t77 = cos(pkin(12));
t81 = sin(qJ(1));
t87 = t81 * t77;
t66 = -t79 * t87 - t86;
t85 = t83 * t77;
t88 = t81 * t73;
t67 = -t79 * t88 + t85;
t68 = t80 * t72 - t82 * t76;
t75 = sin(pkin(6));
t90 = t75 * t81;
t50 = t59 * t90 + t66 * t61 - t67 * t68;
t47 = 0.1e1 / t50 ^ 2;
t58 = t68 * t74;
t60 = t68 * t78;
t48 = -t58 * t90 - t66 * t60 - t67 * t84;
t92 = t48 ^ 2 * t47;
t64 = -t79 * t85 + t88;
t89 = t75 * t83;
t55 = -t64 * t74 + t78 * t89;
t63 = -t75 * t77 * t74 + t79 * t78;
t54 = atan2(t55, t63);
t51 = sin(t54);
t52 = cos(t54);
t45 = t51 * t55 + t52 * t63;
t56 = t66 * t74 - t78 * t90;
t91 = t56 ^ 2 / t45 ^ 2;
t65 = -t79 * t86 - t87;
t62 = 0.1e1 / t63;
t53 = 0.1e1 / (0.1e1 + t55 ^ 2 / t63 ^ 2);
t46 = 0.1e1 / t50;
t42 = 0.1e1 / (0.1e1 + t92);
t1 = [t56 * t62 * t53, 0, 0, 0, 0, 0; (t55 / t45 + (t51 + (t52 * t55 * t62 - t51) * t53) * t91) / (0.1e1 + t91) 0, 0, 0, 0, 0; ((t58 * t89 + t64 * t60 + t65 * t84) * t46 + (t59 * t89 + t64 * t61 - t65 * t68) * t48 * t47) * t42, 0 (t50 * t46 + t92) * t42, 0, 0, 0;];
Ja_rot  = t1;
