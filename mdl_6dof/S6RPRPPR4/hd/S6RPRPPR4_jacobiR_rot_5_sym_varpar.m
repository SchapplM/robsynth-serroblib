% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPPR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:54
% EndTime: 2019-02-26 20:40:54
% DurationCPUTime: 0.07s
% Computational Cost: add. (24->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
t70 = sin(pkin(10));
t72 = sin(qJ(1));
t77 = t72 * t70;
t71 = cos(pkin(10));
t76 = t72 * t71;
t73 = cos(qJ(1));
t75 = t73 * t70;
t74 = t73 * t71;
t69 = pkin(9) + qJ(3);
t68 = cos(t69);
t67 = sin(t69);
t1 = [-t68 * t76 + t75, 0, -t67 * t74, 0, 0, 0; t68 * t74 + t77, 0, -t67 * t76, 0, 0, 0; 0, 0, t68 * t71, 0, 0, 0; -t72 * t67, 0, t73 * t68, 0, 0, 0; t73 * t67, 0, t72 * t68, 0, 0, 0; 0, 0, t67, 0, 0, 0; -t68 * t77 - t74, 0, -t67 * t75, 0, 0, 0; t68 * t75 - t76, 0, -t67 * t77, 0, 0, 0; 0, 0, t68 * t70, 0, 0, 0;];
JR_rot  = t1;
