% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR13
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
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPR13_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:35
% EndTime: 2019-02-26 22:37:36
% DurationCPUTime: 0.03s
% Computational Cost: add. (14->12), mult. (39->18), div. (0->0), fcn. (63->8), ass. (0->18)
t137 = sin(pkin(6));
t141 = sin(qJ(1));
t150 = t141 * t137;
t140 = sin(qJ(2));
t149 = t141 * t140;
t143 = cos(qJ(2));
t148 = t141 * t143;
t144 = cos(qJ(1));
t147 = t144 * t137;
t146 = t144 * t140;
t145 = t144 * t143;
t142 = cos(qJ(3));
t139 = sin(qJ(3));
t138 = cos(pkin(6));
t136 = t137 * t140 * t139 - t138 * t142;
t135 = (-t138 * t149 + t145) * t139 - t142 * t150;
t134 = (t138 * t146 + t148) * t139 + t142 * t147;
t1 = [0, t150, t138 * t148 + t146, t135, 0, -t135; 0, -t147, -t138 * t145 + t149, t134, 0, -t134; 1, t138, -t137 * t143, t136, 0, -t136;];
Jg_rot  = t1;
