% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:52
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPP1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:52:23
% EndTime: 2019-02-22 09:52:23
% DurationCPUTime: 0.07s
% Computational Cost: add. (81->25), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->30)
t145 = qJ(4) + pkin(11);
t143 = sin(t145);
t152 = cos(qJ(3));
t161 = t143 * t152;
t144 = cos(t145);
t160 = t144 * t152;
t147 = sin(pkin(6));
t150 = sin(qJ(3));
t159 = t147 * t150;
t158 = t147 * t152;
t153 = cos(qJ(2));
t157 = t147 * t153;
t149 = cos(pkin(6));
t151 = sin(qJ(2));
t156 = t149 * t151;
t155 = t149 * t153;
t154 = t152 * t153;
t148 = cos(pkin(10));
t146 = sin(pkin(10));
t141 = t149 * t150 + t151 * t158;
t140 = t149 * t152 - t151 * t159;
t139 = -t146 * t156 + t148 * t153;
t138 = t146 * t155 + t148 * t151;
t137 = t146 * t153 + t148 * t156;
t136 = t146 * t151 - t148 * t155;
t135 = t139 * t152 + t146 * t159;
t134 = -t139 * t150 + t146 * t158;
t133 = t137 * t152 - t148 * t159;
t132 = -t137 * t150 - t148 * t158;
t1 = [0, -t138 * t160 + t139 * t143, t134 * t144, -t135 * t143 + t138 * t144, 0, 0; 0, -t136 * t160 + t137 * t143, t132 * t144, -t133 * t143 + t136 * t144, 0, 0; 0 (t143 * t151 + t144 * t154) * t147, t140 * t144, -t141 * t143 - t144 * t157, 0, 0; 0, -t138 * t150, t135, 0, 0, 0; 0, -t136 * t150, t133, 0, 0, 0; 0, t150 * t157, t141, 0, 0, 0; 0, -t138 * t161 - t139 * t144, t134 * t143, t135 * t144 + t138 * t143, 0, 0; 0, -t136 * t161 - t137 * t144, t132 * t143, t133 * t144 + t136 * t143, 0, 0; 0 (t143 * t154 - t144 * t151) * t147, t140 * t143, t141 * t144 - t143 * t157, 0, 0;];
JR_rot  = t1;
