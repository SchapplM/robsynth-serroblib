% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPPRR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:50
% EndTime: 2019-02-26 19:45:50
% DurationCPUTime: 0.05s
% Computational Cost: add. (44->22), mult. (122->44), div. (0->0), fcn. (178->10), ass. (0->25)
t101 = sin(qJ(5));
t97 = sin(pkin(6));
t108 = t101 * t97;
t103 = cos(qJ(5));
t107 = t103 * t97;
t100 = cos(pkin(6));
t102 = sin(qJ(2));
t106 = t100 * t102;
t104 = cos(qJ(2));
t105 = t100 * t104;
t96 = sin(pkin(10));
t99 = cos(pkin(10));
t90 = t96 * t102 - t99 * t105;
t91 = t96 * t104 + t99 * t106;
t95 = sin(pkin(11));
t98 = cos(pkin(11));
t85 = t90 * t95 + t91 * t98;
t92 = t99 * t102 + t96 * t105;
t93 = t99 * t104 - t96 * t106;
t87 = t92 * t95 + t93 * t98;
t89 = (t102 * t98 - t104 * t95) * t97;
t88 = (t102 * t95 + t104 * t98) * t97;
t86 = -t92 * t98 + t93 * t95;
t84 = -t90 * t98 + t91 * t95;
t1 = [0, t86 * t103, 0, 0, -t87 * t101 - t96 * t107, 0; 0, t84 * t103, 0, 0, -t85 * t101 + t99 * t107, 0; 0, t88 * t103, 0, 0, -t100 * t103 - t89 * t101, 0; 0, -t86 * t101, 0, 0, -t87 * t103 + t96 * t108, 0; 0, -t84 * t101, 0, 0, -t85 * t103 - t99 * t108, 0; 0, -t88 * t101, 0, 0, t100 * t101 - t89 * t103, 0; 0, -t87, 0, 0, 0, 0; 0, -t85, 0, 0, 0, 0; 0, -t89, 0, 0, 0, 0;];
JR_rot  = t1;
