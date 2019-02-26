% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPPRR1_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobiR_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:35
% DurationCPUTime: 0.10s
% Computational Cost: add. (24->9), mult. (66->22), div. (0->0), fcn. (96->10), ass. (0->18)
t74 = sin(pkin(11));
t78 = cos(pkin(11));
t81 = sin(qJ(2));
t82 = cos(qJ(2));
t83 = t74 * t82 + t78 * t81;
t71 = t74 * t81 - t78 * t82;
t80 = cos(pkin(6));
t79 = cos(pkin(10));
t77 = cos(pkin(12));
t76 = sin(pkin(6));
t75 = sin(pkin(10));
t73 = sin(pkin(12));
t70 = t83 * t80;
t69 = t71 * t80;
t68 = t71 * t76;
t67 = t69 * t75 - t79 * t83;
t66 = -t69 * t79 - t75 * t83;
t1 = [0, t67 * t77, 0, 0, 0, 0; 0, t66 * t77, 0, 0, 0, 0; 0, -t68 * t77, 0, 0, 0, 0; 0, -t67 * t73, 0, 0, 0, 0; 0, -t66 * t73, 0, 0, 0, 0; 0, t68 * t73, 0, 0, 0, 0; 0, -t70 * t75 - t71 * t79, 0, 0, 0, 0; 0, t70 * t79 - t71 * t75, 0, 0, 0, 0; 0, t83 * t76, 0, 0, 0, 0;];
JR_rot  = t1;
