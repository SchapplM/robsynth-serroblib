% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR2_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:12
% EndTime: 2019-02-26 19:54:12
% DurationCPUTime: 0.04s
% Computational Cost: add. (28->10), mult. (78->24), div. (0->0), fcn. (117->10), ass. (0->20)
t130 = sin(pkin(11));
t131 = sin(pkin(6));
t142 = t130 * t131;
t133 = cos(pkin(11));
t141 = t133 * t131;
t129 = sin(pkin(12));
t132 = cos(pkin(12));
t136 = sin(qJ(2));
t138 = cos(qJ(2));
t140 = t138 * t129 + t136 * t132;
t139 = t136 * t129 - t138 * t132;
t137 = cos(qJ(4));
t135 = sin(qJ(4));
t134 = cos(pkin(6));
t126 = t140 * t134;
t125 = t139 * t134;
t124 = t140 * t135 * t131 - t134 * t137;
t123 = (-t130 * t126 - t133 * t139) * t135 - t137 * t142;
t122 = (t133 * t126 - t130 * t139) * t135 + t137 * t141;
t1 = [0, t142, 0, -t130 * t125 + t133 * t140, t123, t123; 0, -t141, 0, t133 * t125 + t130 * t140, t122, t122; 0, t134, 0, t139 * t131, t124, t124;];
Jg_rot  = t1;
