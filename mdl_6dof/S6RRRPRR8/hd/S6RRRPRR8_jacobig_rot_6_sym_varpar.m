% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR8_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:57
% EndTime: 2019-02-26 22:19:57
% DurationCPUTime: 0.03s
% Computational Cost: add. (26->10), mult. (39->18), div. (0->0), fcn. (63->8), ass. (0->19)
t128 = sin(pkin(6));
t131 = sin(qJ(1));
t139 = t131 * t128;
t130 = sin(qJ(2));
t138 = t131 * t130;
t132 = cos(qJ(2));
t137 = t131 * t132;
t133 = cos(qJ(1));
t136 = t133 * t128;
t135 = t133 * t130;
t134 = t133 * t132;
t129 = cos(pkin(6));
t127 = qJ(3) + pkin(12);
t126 = cos(t127);
t125 = sin(t127);
t124 = t128 * t130 * t125 - t129 * t126;
t123 = (-t129 * t138 + t134) * t125 - t126 * t139;
t122 = (t129 * t135 + t137) * t125 + t126 * t136;
t1 = [0, t139, t129 * t137 + t135, 0, t123, t123; 0, -t136, -t129 * t134 + t138, 0, t122, t122; 1, t129, -t128 * t132, 0, t124, t124;];
Jg_rot  = t1;
