% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPP4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:03
% EndTime: 2019-02-26 20:58:03
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t96 = qJ(4) + pkin(10);
t92 = sin(t96);
t97 = sin(qJ(1));
t104 = t97 * t92;
t95 = pkin(9) + qJ(3);
t93 = cos(t95);
t103 = t97 * t93;
t94 = cos(t96);
t102 = t97 * t94;
t98 = cos(qJ(1));
t101 = t98 * t92;
t100 = t98 * t93;
t99 = t98 * t94;
t91 = sin(t95);
t90 = t93 * t99 + t104;
t89 = t92 * t100 - t102;
t88 = t93 * t102 - t101;
t87 = -t92 * t103 - t99;
t1 = [-t88, 0, -t91 * t99, -t89, 0, 0; t90, 0, -t91 * t102, t87, 0, 0; 0, 0, t93 * t94, -t91 * t92, 0, 0; -t97 * t91, 0, t100, 0, 0, 0; t98 * t91, 0, t103, 0, 0, 0; 0, 0, t91, 0, 0, 0; t87, 0, -t91 * t101, t90, 0, 0; t89, 0, -t91 * t104, t88, 0, 0; 0, 0, t93 * t92, t91 * t94, 0, 0;];
JR_rot  = t1;
