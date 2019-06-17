% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-12 14:37
% Revision: aab8d7cd0cba739f5e0ec8d53b8419901d1154b0 (2019-06-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRRR1_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiR_rot_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiR_rot_4_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-12 14:37:34
% EndTime: 2019-06-12 14:37:34
% DurationCPUTime: 0.05s
% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
t24 = sin(qJ(3));
t25 = sin(qJ(1));
t35 = t25 * t24;
t26 = cos(qJ(4));
t34 = t25 * t26;
t23 = sin(qJ(4));
t27 = cos(qJ(3));
t33 = t27 * t23;
t32 = t27 * t26;
t28 = cos(qJ(1));
t31 = t28 * t24;
t30 = t28 * t26;
t29 = t28 * t27;
t22 = t25 * t23 + t26 * t29;
t21 = -t23 * t29 + t34;
t20 = t28 * t23 - t25 * t32;
t19 = t25 * t33 + t30;
t1 = [t20, 0, -t24 * t30, t21, 0; t22, 0, -t24 * t34, -t19, 0; 0, 0, t32, -t24 * t23, 0; t19, 0, t23 * t31, -t22, 0; t21, 0, t23 * t35, t20, 0; 0, 0, -t33, -t24 * t26, 0; -t35, 0, t29, 0, 0; t31, 0, t25 * t27, 0, 0; 0, 0, t24, 0, 0;];
JR_rot  = t1;
