% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:45
% EndTime: 2019-02-26 19:48:45
% DurationCPUTime: 0.04s
% Computational Cost: add. (37->14), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->20)
t91 = sin(pkin(10));
t92 = sin(pkin(6));
t102 = t91 * t92;
t93 = cos(pkin(10));
t101 = t92 * t93;
t95 = sin(qJ(2));
t100 = t92 * t95;
t96 = cos(qJ(2));
t99 = t92 * t96;
t94 = cos(pkin(6));
t98 = t94 * t95;
t97 = t94 * t96;
t90 = pkin(11) + qJ(4);
t89 = cos(t90);
t88 = sin(t90);
t87 = -t91 * t98 + t93 * t96;
t86 = -t91 * t97 - t93 * t95;
t85 = t91 * t96 + t93 * t98;
t84 = -t91 * t95 + t93 * t97;
t1 = [0, t87, 0, 0, 0, 0; 0, t85, 0, 0, 0, 0; 0, t100, 0, 0, 0, 0; 0, -t86 * t89, 0, -t89 * t102 + t87 * t88, 0, 0; 0, -t84 * t89, 0, t89 * t101 + t85 * t88, 0, 0; 0, -t89 * t99, 0, t88 * t100 - t94 * t89, 0, 0; 0, t86 * t88, 0, t88 * t102 + t87 * t89, 0, 0; 0, t84 * t88, 0, -t88 * t101 + t85 * t89, 0, 0; 0, t88 * t99, 0, t89 * t100 + t94 * t88, 0, 0;];
JR_rot  = t1;
