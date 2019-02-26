% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPR10_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:53
% EndTime: 2019-02-26 22:35:53
% DurationCPUTime: 0.03s
% Computational Cost: add. (18->11), mult. (31->18), div. (0->0), fcn. (52->8), ass. (0->19)
t129 = sin(pkin(6));
t133 = cos(qJ(2));
t141 = t129 * t133;
t132 = sin(qJ(1));
t140 = t132 * t129;
t131 = sin(qJ(2));
t139 = t132 * t131;
t138 = t132 * t133;
t134 = cos(qJ(1));
t137 = t134 * t129;
t136 = t134 * t131;
t135 = t134 * t133;
t130 = cos(pkin(6));
t128 = qJ(3) + qJ(4);
t127 = cos(t128);
t126 = sin(t128);
t125 = t130 * t138 + t136;
t124 = -t130 * t135 + t139;
t1 = [0, t140, t125, t125, 0 (-t130 * t139 + t135) * t127 + t126 * t140; 0, -t137, t124, t124, 0 (t130 * t136 + t138) * t127 - t126 * t137; 1, t130, -t141, -t141, 0, t129 * t131 * t127 + t130 * t126;];
Jg_rot  = t1;
