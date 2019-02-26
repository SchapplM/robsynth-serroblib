% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:46
% EndTime: 2019-02-26 19:58:46
% DurationCPUTime: 0.04s
% Computational Cost: add. (37->14), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->20)
t94 = sin(pkin(10));
t95 = sin(pkin(6));
t105 = t94 * t95;
t96 = cos(pkin(10));
t104 = t95 * t96;
t98 = sin(qJ(2));
t103 = t95 * t98;
t99 = cos(qJ(2));
t102 = t95 * t99;
t97 = cos(pkin(6));
t101 = t97 * t98;
t100 = t97 * t99;
t93 = qJ(3) + pkin(11);
t92 = cos(t93);
t91 = sin(t93);
t90 = -t94 * t101 + t96 * t99;
t89 = -t94 * t100 - t96 * t98;
t88 = t96 * t101 + t94 * t99;
t87 = t96 * t100 - t94 * t98;
t1 = [0, t90, 0, 0, 0, 0; 0, t88, 0, 0, 0, 0; 0, t103, 0, 0, 0, 0; 0, -t89 * t92, -t92 * t105 + t90 * t91, 0, 0, 0; 0, -t87 * t92, t92 * t104 + t88 * t91, 0, 0, 0; 0, -t92 * t102, t91 * t103 - t97 * t92, 0, 0, 0; 0, t89 * t91, t91 * t105 + t90 * t92, 0, 0, 0; 0, t87 * t91, -t91 * t104 + t88 * t92, 0, 0, 0; 0, t91 * t102, t92 * t103 + t97 * t91, 0, 0, 0;];
JR_rot  = t1;
