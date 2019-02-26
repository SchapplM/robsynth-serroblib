% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:47
% EndTime: 2019-02-26 21:02:47
% DurationCPUTime: 0.03s
% Computational Cost: add. (59->12), mult. (38->18), div. (0->0), fcn. (66->6), ass. (0->20)
t88 = pkin(10) + qJ(3) + qJ(4);
t87 = cos(t88);
t89 = sin(pkin(11));
t99 = t87 * t89;
t91 = sin(qJ(1));
t98 = t91 * t89;
t90 = cos(pkin(11));
t97 = t91 * t90;
t92 = cos(qJ(1));
t96 = t92 * t89;
t95 = t92 * t90;
t86 = sin(t88);
t94 = t86 * t97;
t93 = t86 * t95;
t85 = t92 * t87;
t84 = t91 * t87;
t83 = t87 * t90;
t82 = t86 * t96;
t81 = t86 * t98;
t1 = [-t87 * t97 + t96, 0, -t93, -t93, 0, 0; t87 * t95 + t98, 0, -t94, -t94, 0, 0; 0, 0, t83, t83, 0, 0; t87 * t98 + t95, 0, t82, t82, 0, 0; -t87 * t96 + t97, 0, t81, t81, 0, 0; 0, 0, -t99, -t99, 0, 0; -t91 * t86, 0, t85, t85, 0, 0; t92 * t86, 0, t84, t84, 0, 0; 0, 0, t86, t86, 0, 0;];
JR_rot  = t1;
