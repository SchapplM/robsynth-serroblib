% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRP6_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:50
% EndTime: 2019-02-26 21:48:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
t160 = sin(pkin(6));
t165 = sin(qJ(1));
t172 = t165 * t160;
t168 = cos(qJ(1));
t171 = t168 * t160;
t159 = sin(pkin(11));
t161 = cos(pkin(11));
t164 = sin(qJ(2));
t167 = cos(qJ(2));
t170 = t167 * t159 + t164 * t161;
t169 = t164 * t159 - t167 * t161;
t166 = cos(qJ(4));
t163 = sin(qJ(4));
t162 = cos(pkin(6));
t156 = t170 * t162;
t155 = t169 * t162;
t1 = [0, t172, 0, -t165 * t155 + t168 * t170 (-t165 * t156 - t168 * t169) * t163 - t166 * t172, 0; 0, -t171, 0, t168 * t155 + t165 * t170 (t168 * t156 - t165 * t169) * t163 + t166 * t171, 0; 1, t162, 0, t169 * t160, t170 * t163 * t160 - t162 * t166, 0;];
Jg_rot  = t1;
