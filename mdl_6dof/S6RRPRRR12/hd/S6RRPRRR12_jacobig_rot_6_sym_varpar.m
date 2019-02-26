% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR12_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:00:24
% EndTime: 2019-02-26 22:00:24
% DurationCPUTime: 0.03s
% Computational Cost: add. (16->9), mult. (31->18), div. (0->0), fcn. (52->8), ass. (0->19)
t129 = sin(pkin(6));
t132 = sin(qJ(1));
t140 = t132 * t129;
t131 = sin(qJ(2));
t139 = t132 * t131;
t133 = cos(qJ(2));
t138 = t132 * t133;
t134 = cos(qJ(1));
t137 = t134 * t129;
t136 = t134 * t131;
t135 = t134 * t133;
t130 = cos(pkin(6));
t128 = qJ(4) + qJ(5);
t127 = cos(t128);
t126 = sin(t128);
t125 = t129 * t131;
t124 = -t130 * t139 + t135;
t123 = t130 * t136 + t138;
t1 = [0, t140, 0, t124, t124, t126 * t140 - (t130 * t138 + t136) * t127; 0, -t137, 0, t123, t123, -t126 * t137 - (-t130 * t135 + t139) * t127; 1, t130, 0, t125, t125, t129 * t133 * t127 + t130 * t126;];
Jg_rot  = t1;
