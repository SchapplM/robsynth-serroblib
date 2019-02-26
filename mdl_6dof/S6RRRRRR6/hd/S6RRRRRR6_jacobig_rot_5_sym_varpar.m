% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR6_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:04
% EndTime: 2019-02-26 22:50:04
% DurationCPUTime: 0.03s
% Computational Cost: add. (18->11), mult. (31->18), div. (0->0), fcn. (52->8), ass. (0->19)
t130 = sin(pkin(6));
t134 = cos(qJ(2));
t142 = t130 * t134;
t133 = sin(qJ(1));
t141 = t133 * t130;
t132 = sin(qJ(2));
t140 = t133 * t132;
t139 = t133 * t134;
t135 = cos(qJ(1));
t138 = t135 * t130;
t137 = t135 * t132;
t136 = t135 * t134;
t131 = cos(pkin(6));
t129 = qJ(3) + qJ(4);
t128 = cos(t129);
t127 = sin(t129);
t126 = t131 * t139 + t137;
t125 = -t131 * t136 + t140;
t1 = [0, t141, t126, t126 (-t131 * t140 + t136) * t127 - t128 * t141, 0; 0, -t138, t125, t125 (t131 * t137 + t139) * t127 + t128 * t138, 0; 1, t131, -t142, -t142, t130 * t132 * t127 - t131 * t128, 0;];
Jg_rot  = t1;
