% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:29
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:29:20
% EndTime: 2019-02-22 09:29:20
% DurationCPUTime: 0.10s
% Computational Cost: add. (153->29), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->34)
t156 = sin(pkin(11));
t159 = cos(pkin(11));
t163 = sin(qJ(2));
t165 = cos(qJ(2));
t149 = t163 * t156 - t165 * t159;
t155 = qJ(4) + pkin(12);
t154 = cos(t155);
t162 = sin(qJ(6));
t175 = t154 * t162;
t164 = cos(qJ(6));
t174 = t154 * t164;
t157 = sin(pkin(10));
t158 = sin(pkin(6));
t173 = t157 * t158;
t160 = cos(pkin(10));
t172 = t158 * t160;
t161 = cos(pkin(6));
t167 = t165 * t156 + t163 * t159;
t148 = t167 * t161;
t169 = t160 * t148 - t157 * t149;
t168 = -t157 * t148 - t160 * t149;
t166 = t149 * t161;
t153 = sin(t155);
t147 = t167 * t158;
t146 = t149 * t158;
t144 = t147 * t154 + t161 * t153;
t143 = -t147 * t153 + t161 * t154;
t141 = t157 * t166 - t160 * t167;
t138 = -t157 * t167 - t160 * t166;
t136 = t153 * t173 + t154 * t168;
t135 = -t153 * t168 + t154 * t173;
t134 = -t153 * t172 + t154 * t169;
t133 = -t153 * t169 - t154 * t172;
t1 = [0, t141 * t174 + t162 * t168, 0, t135 * t164, 0, -t136 * t162 - t141 * t164; 0, t138 * t174 + t162 * t169, 0, t133 * t164, 0, -t134 * t162 - t138 * t164; 0, -t146 * t174 + t147 * t162, 0, t143 * t164, 0, -t144 * t162 + t146 * t164; 0, -t141 * t175 + t164 * t168, 0, -t135 * t162, 0, -t136 * t164 + t141 * t162; 0, -t138 * t175 + t164 * t169, 0, -t133 * t162, 0, -t134 * t164 + t138 * t162; 0, t146 * t175 + t147 * t164, 0, -t143 * t162, 0, -t144 * t164 - t146 * t162; 0, t141 * t153, 0, t136, 0, 0; 0, t138 * t153, 0, t134, 0, 0; 0, -t146 * t153, 0, t144, 0, 0;];
JR_rot  = t1;
