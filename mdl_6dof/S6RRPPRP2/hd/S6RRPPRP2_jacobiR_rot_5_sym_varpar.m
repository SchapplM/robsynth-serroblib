% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRP2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:25:21
% EndTime: 2019-02-26 21:25:22
% DurationCPUTime: 0.03s
% Computational Cost: add. (33->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
t75 = sin(qJ(5));
t76 = sin(qJ(1));
t82 = t76 * t75;
t77 = cos(qJ(5));
t81 = t76 * t77;
t78 = cos(qJ(1));
t80 = t78 * t75;
t79 = t78 * t77;
t74 = qJ(2) + pkin(9);
t73 = cos(t74);
t72 = sin(t74);
t71 = -t72 * t82 + t79;
t70 = t72 * t81 + t80;
t69 = t72 * t80 + t81;
t68 = t72 * t79 - t82;
t1 = [t71, t73 * t80, 0, 0, t68, 0; t69, t73 * t82, 0, 0, t70, 0; 0, t72 * t75, 0, 0, -t73 * t77, 0; -t70, t73 * t79, 0, 0, -t69, 0; t68, t73 * t81, 0, 0, t71, 0; 0, t72 * t77, 0, 0, t73 * t75, 0; -t76 * t73, -t78 * t72, 0, 0, 0, 0; t78 * t73, -t76 * t72, 0, 0, 0, 0; 0, t73, 0, 0, 0, 0;];
JR_rot  = t1;
