% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:06
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPPRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:06:01
% EndTime: 2019-02-22 10:06:01
% DurationCPUTime: 0.02s
% Computational Cost: add. (18->9), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
t14 = qJ(1) + pkin(9);
t12 = sin(t14);
t15 = sin(qJ(5));
t19 = t12 * t15;
t16 = cos(qJ(5));
t18 = t12 * t16;
t13 = cos(t14);
t17 = t13 * t15;
t11 = t13 * t16;
t1 = [-t19, 0, 0, 0, t11, 0; t17, 0, 0, 0, t18, 0; 0, 0, 0, 0, -t15, 0; -t18, 0, 0, 0, -t17, 0; t11, 0, 0, 0, -t19, 0; 0, 0, 0, 0, -t16, 0; -t13, 0, 0, 0, 0, 0; -t12, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;