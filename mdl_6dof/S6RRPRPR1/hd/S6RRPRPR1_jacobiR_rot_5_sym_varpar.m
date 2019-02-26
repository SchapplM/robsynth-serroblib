% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:40
% EndTime: 2019-02-26 21:37:40
% DurationCPUTime: 0.06s
% Computational Cost: add. (59->12), mult. (38->18), div. (0->0), fcn. (66->6), ass. (0->20)
t90 = qJ(2) + pkin(10) + qJ(4);
t89 = cos(t90);
t91 = sin(pkin(11));
t101 = t89 * t91;
t93 = sin(qJ(1));
t100 = t93 * t91;
t92 = cos(pkin(11));
t99 = t93 * t92;
t94 = cos(qJ(1));
t98 = t94 * t91;
t97 = t94 * t92;
t88 = sin(t90);
t96 = t88 * t99;
t95 = t88 * t97;
t87 = t94 * t89;
t86 = t93 * t89;
t85 = t89 * t92;
t84 = t88 * t98;
t83 = t88 * t100;
t1 = [-t89 * t99 + t98, -t95, 0, -t95, 0, 0; t89 * t97 + t100, -t96, 0, -t96, 0, 0; 0, t85, 0, t85, 0, 0; t89 * t100 + t97, t84, 0, t84, 0, 0; -t89 * t98 + t99, t83, 0, t83, 0, 0; 0, -t101, 0, -t101, 0, 0; -t93 * t88, t87, 0, t87, 0, 0; t94 * t88, t86, 0, t86, 0, 0; 0, t88, 0, t88, 0, 0;];
JR_rot  = t1;
