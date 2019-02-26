% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:47
% EndTime: 2019-02-26 19:49:47
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->18), mult. (57->31), div. (0->0), fcn. (88->8), ass. (0->18)
t82 = sin(pkin(6));
t85 = sin(qJ(4));
t93 = t82 * t85;
t87 = cos(qJ(4));
t92 = t82 * t87;
t88 = cos(qJ(2));
t91 = t82 * t88;
t84 = cos(pkin(6));
t86 = sin(qJ(2));
t90 = t84 * t86;
t89 = t84 * t88;
t83 = cos(pkin(10));
t81 = sin(pkin(10));
t80 = -t81 * t90 + t83 * t88;
t79 = t81 * t89 + t83 * t86;
t78 = t81 * t88 + t83 * t90;
t77 = t81 * t86 - t83 * t89;
t1 = [0, -t79, 0, 0, 0, 0; 0, -t77, 0, 0, 0, 0; 0, t91, 0, 0, 0, 0; 0, -t80 * t85, 0, -t79 * t87 + t81 * t93, 0, 0; 0, -t78 * t85, 0, -t77 * t87 - t83 * t93, 0, 0; 0, -t86 * t93, 0, t84 * t85 + t87 * t91, 0, 0; 0, -t80 * t87, 0, t79 * t85 + t81 * t92, 0, 0; 0, -t78 * t87, 0, t77 * t85 - t83 * t92, 0, 0; 0, -t86 * t92, 0, t84 * t87 - t85 * t91, 0, 0;];
JR_rot  = t1;
