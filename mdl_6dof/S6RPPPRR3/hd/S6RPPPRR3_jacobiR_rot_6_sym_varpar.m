% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPPRR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:23:47
% EndTime: 2019-02-26 20:23:47
% DurationCPUTime: 0.04s
% Computational Cost: add. (57->16), mult. (88->24), div. (0->0), fcn. (141->8), ass. (0->20)
t102 = cos(qJ(1));
t101 = sin(qJ(1));
t90 = pkin(10) + qJ(5);
t88 = sin(t90);
t91 = sin(qJ(6));
t100 = t88 * t91;
t92 = cos(qJ(6));
t99 = t88 * t92;
t89 = cos(t90);
t98 = t89 * t91;
t97 = t89 * t92;
t96 = cos(pkin(9));
t95 = sin(pkin(9));
t83 = -t101 * t95 - t102 * t96;
t84 = -t101 * t96 + t102 * t95;
t94 = t83 * t91 + t84 * t97;
t93 = -t83 * t92 + t84 * t98;
t82 = -t83 * t97 + t84 * t91;
t81 = t83 * t98 + t84 * t92;
t1 = [t94, 0, 0, 0, t83 * t99, t81; t82, 0, 0, 0, t84 * t99, t93; 0, 0, 0, 0, -t97, t100; -t93, 0, 0, 0, -t83 * t100, -t82; t81, 0, 0, 0, -t84 * t100, t94; 0, 0, 0, 0, t98, t99; t84 * t88, 0, 0, 0, -t83 * t89, 0; -t83 * t88, 0, 0, 0, -t84 * t89, 0; 0, 0, 0, 0, -t88, 0;];
JR_rot  = t1;
