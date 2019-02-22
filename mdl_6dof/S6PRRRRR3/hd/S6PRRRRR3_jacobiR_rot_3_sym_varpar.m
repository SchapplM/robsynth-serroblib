% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:03
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR3_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_jacobiR_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_jacobiR_rot_3_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:03:02
% EndTime: 2019-02-22 10:03:02
% DurationCPUTime: 0.06s
% Computational Cost: add. (19->13), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
t65 = sin(pkin(6));
t68 = sin(qJ(3));
t77 = t65 * t68;
t69 = sin(qJ(2));
t76 = t65 * t69;
t70 = cos(qJ(3));
t75 = t65 * t70;
t71 = cos(qJ(2));
t74 = t65 * t71;
t67 = cos(pkin(6));
t73 = t67 * t69;
t72 = t67 * t71;
t66 = cos(pkin(12));
t64 = sin(pkin(12));
t63 = -t64 * t73 + t66 * t71;
t62 = -t64 * t72 - t66 * t69;
t61 = t64 * t71 + t66 * t73;
t60 = -t64 * t69 + t66 * t72;
t1 = [0, t62 * t70, -t63 * t68 + t64 * t75, 0, 0, 0; 0, t60 * t70, -t61 * t68 - t66 * t75, 0, 0, 0; 0, t70 * t74, t67 * t70 - t68 * t76, 0, 0, 0; 0, -t62 * t68, -t63 * t70 - t64 * t77, 0, 0, 0; 0, -t60 * t68, -t61 * t70 + t66 * t77, 0, 0, 0; 0, -t68 * t74, -t67 * t68 - t69 * t75, 0, 0, 0; 0, t63, 0, 0, 0, 0; 0, t61, 0, 0, 0, 0; 0, t76, 0, 0, 0, 0;];
JR_rot  = t1;
