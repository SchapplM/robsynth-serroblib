% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRPR9_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:29
% EndTime: 2019-02-26 21:05:29
% DurationCPUTime: 0.04s
% Computational Cost: add. (15->13), mult. (46->31), div. (0->0), fcn. (67->10), ass. (0->19)
t131 = sin(pkin(6));
t136 = sin(qJ(1));
t144 = t131 * t136;
t138 = cos(qJ(1));
t143 = t131 * t138;
t129 = sin(pkin(12));
t142 = t136 * t129;
t132 = cos(pkin(12));
t141 = t136 * t132;
t140 = t138 * t129;
t139 = t138 * t132;
t137 = cos(qJ(3));
t135 = sin(qJ(3));
t134 = cos(pkin(6));
t133 = cos(pkin(7));
t130 = sin(pkin(7));
t128 = -t134 * t141 - t140;
t127 = t134 * t139 - t142;
t1 = [0, 0, -t128 * t130 + t133 * t144 (-t134 * t142 + t139) * t135 + (-t128 * t133 - t130 * t144) * t137, 0, 0; 0, 0, -t127 * t130 - t133 * t143 (t134 * t140 + t141) * t135 + (-t127 * t133 + t130 * t143) * t137, 0, 0; 1, 0, -t131 * t132 * t130 + t134 * t133, -t134 * t130 * t137 + (-t132 * t133 * t137 + t129 * t135) * t131, 0, 0;];
Jg_rot  = t1;
