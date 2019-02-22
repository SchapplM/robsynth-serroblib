% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:41
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR7_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiR_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiR_rot_3_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:41:40
% EndTime: 2019-02-22 09:41:40
% DurationCPUTime: 0.06s
% Computational Cost: add. (18->12), mult. (56->31), div. (0->0), fcn. (81->10), ass. (0->20)
t64 = sin(pkin(14));
t70 = cos(pkin(7));
t78 = t64 * t70;
t68 = cos(pkin(14));
t77 = t68 * t70;
t72 = sin(qJ(2));
t76 = t70 * t72;
t71 = cos(pkin(6));
t75 = t71 * t72;
t73 = cos(qJ(2));
t74 = t71 * t73;
t69 = cos(pkin(13));
t67 = sin(pkin(6));
t66 = sin(pkin(7));
t65 = sin(pkin(13));
t63 = t65 * t75 - t69 * t73;
t62 = -t65 * t74 - t69 * t72;
t61 = -t65 * t73 - t69 * t75;
t60 = -t65 * t72 + t69 * t74;
t1 = [0, t62 * t68 + t63 * t78, 0, 0, 0, 0; 0, t60 * t68 + t61 * t78, 0, 0, 0, 0; 0 (-t64 * t76 + t68 * t73) * t67, 0, 0, 0, 0; 0, -t62 * t64 + t63 * t77, 0, 0, 0, 0; 0, -t60 * t64 + t61 * t77, 0, 0, 0, 0; 0 (-t64 * t73 - t68 * t76) * t67, 0, 0, 0, 0; 0, -t63 * t66, 0, 0, 0, 0; 0, -t61 * t66, 0, 0, 0, 0; 0, t67 * t72 * t66, 0, 0, 0, 0;];
JR_rot  = t1;
