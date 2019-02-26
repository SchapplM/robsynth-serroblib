% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRPR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:16
% EndTime: 2019-02-26 20:26:16
% DurationCPUTime: 0.03s
% Computational Cost: add. (57->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t77 = qJ(1) + pkin(9);
t73 = sin(t77);
t78 = sin(qJ(6));
t85 = t73 * t78;
t79 = cos(qJ(6));
t84 = t73 * t79;
t76 = pkin(10) + qJ(4);
t74 = cos(t76);
t83 = t74 * t78;
t82 = t74 * t79;
t75 = cos(t77);
t81 = t75 * t78;
t80 = t75 * t79;
t72 = sin(t76);
t71 = -t72 * t85 + t80;
t70 = t72 * t84 + t81;
t69 = t72 * t81 + t84;
t68 = t72 * t80 - t85;
t1 = [t71, 0, 0, t74 * t81, 0, t68; t69, 0, 0, t73 * t83, 0, t70; 0, 0, 0, t72 * t78, 0, -t82; -t70, 0, 0, t74 * t80, 0, -t69; t68, 0, 0, t73 * t82, 0, t71; 0, 0, 0, t72 * t79, 0, t83; -t73 * t74, 0, 0, -t75 * t72, 0, 0; t75 * t74, 0, 0, -t73 * t72, 0, 0; 0, 0, 0, t74, 0, 0;];
JR_rot  = t1;
