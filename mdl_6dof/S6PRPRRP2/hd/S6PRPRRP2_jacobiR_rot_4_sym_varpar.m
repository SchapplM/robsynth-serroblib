% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRP2_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobiR_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:50:54
% EndTime: 2019-02-26 19:50:54
% DurationCPUTime: 0.05s
% Computational Cost: add. (44->15), mult. (122->36), div. (0->0), fcn. (178->10), ass. (0->23)
t93 = sin(pkin(6));
t97 = sin(qJ(4));
t104 = t93 * t97;
t99 = cos(qJ(4));
t103 = t93 * t99;
t100 = cos(qJ(2));
t91 = sin(pkin(11));
t94 = cos(pkin(11));
t98 = sin(qJ(2));
t102 = t100 * t94 - t98 * t91;
t101 = t100 * t91 + t98 * t94;
t96 = cos(pkin(6));
t87 = t101 * t96;
t92 = sin(pkin(10));
t95 = cos(pkin(10));
t81 = t102 * t92 + t95 * t87;
t83 = t102 * t95 - t92 * t87;
t86 = t102 * t96;
t85 = t101 * t93;
t84 = t102 * t93;
t82 = -t101 * t95 - t92 * t86;
t80 = -t101 * t92 + t95 * t86;
t1 = [0, t82 * t99, 0, t92 * t103 - t83 * t97, 0, 0; 0, t80 * t99, 0, -t95 * t103 - t81 * t97, 0, 0; 0, t84 * t99, 0, -t85 * t97 + t96 * t99, 0, 0; 0, -t82 * t97, 0, -t92 * t104 - t83 * t99, 0, 0; 0, -t80 * t97, 0, t95 * t104 - t81 * t99, 0, 0; 0, -t84 * t97, 0, -t85 * t99 - t96 * t97, 0, 0; 0, t83, 0, 0, 0, 0; 0, t81, 0, 0, 0, 0; 0, t85, 0, 0, 0, 0;];
JR_rot  = t1;
