% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPPRR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:24:53
% EndTime: 2019-02-26 20:24:53
% DurationCPUTime: 0.06s
% Computational Cost: add. (37->14), mult. (88->24), div. (0->0), fcn. (141->8), ass. (0->19)
t89 = sin(qJ(1));
t77 = sin(qJ(6));
t78 = sin(qJ(5));
t88 = t78 * t77;
t79 = cos(qJ(6));
t87 = t78 * t79;
t80 = cos(qJ(5));
t86 = t80 * t77;
t85 = t80 * t79;
t84 = sin(pkin(9));
t76 = cos(pkin(9));
t81 = cos(qJ(1));
t73 = t81 * t76 - t89 * t84;
t74 = t89 * t76 + t81 * t84;
t83 = t73 * t85 + t74 * t77;
t82 = t73 * t86 - t74 * t79;
t72 = -t73 * t77 + t74 * t85;
t71 = -t73 * t79 - t74 * t86;
t1 = [t83, 0, 0, 0, -t74 * t87, t71; t72, 0, 0, 0, t73 * t87, t82; 0, 0, 0, 0, t85, -t88; -t82, 0, 0, 0, t74 * t88, -t72; t71, 0, 0, 0, -t73 * t88, t83; 0, 0, 0, 0, -t86, -t87; t73 * t78, 0, 0, 0, t74 * t80, 0; t74 * t78, 0, 0, 0, -t73 * t80, 0; 0, 0, 0, 0, t78, 0;];
JR_rot  = t1;
