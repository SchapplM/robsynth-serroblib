% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:59:23
% EndTime: 2019-02-26 19:59:23
% DurationCPUTime: 0.06s
% Computational Cost: add. (20->16), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
t76 = sin(pkin(6));
t79 = sin(qJ(3));
t88 = t76 * t79;
t80 = sin(qJ(2));
t87 = t76 * t80;
t81 = cos(qJ(3));
t86 = t76 * t81;
t82 = cos(qJ(2));
t85 = t76 * t82;
t78 = cos(pkin(6));
t84 = t78 * t80;
t83 = t78 * t82;
t77 = cos(pkin(10));
t75 = sin(pkin(10));
t74 = -t75 * t84 + t77 * t82;
t73 = -t75 * t83 - t77 * t80;
t72 = t75 * t82 + t77 * t84;
t71 = -t75 * t80 + t77 * t83;
t1 = [0, t73 * t79, t74 * t81 + t75 * t88, 0, 0, 0; 0, t71 * t79, t72 * t81 - t77 * t88, 0, 0, 0; 0, t79 * t85, t78 * t79 + t80 * t86, 0, 0, 0; 0, -t73 * t81, t74 * t79 - t75 * t86, 0, 0, 0; 0, -t71 * t81, t72 * t79 + t77 * t86, 0, 0, 0; 0, -t81 * t85, -t78 * t81 + t79 * t87, 0, 0, 0; 0, -t74, 0, 0, 0, 0; 0, -t72, 0, 0, 0, 0; 0, -t87, 0, 0, 0, 0;];
JR_rot  = t1;
