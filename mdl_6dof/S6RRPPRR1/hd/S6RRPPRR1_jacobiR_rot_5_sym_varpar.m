% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:12
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:12:46
% EndTime: 2019-02-22 11:12:46
% DurationCPUTime: 0.04s
% Computational Cost: add. (50->13), mult. (56->8), div. (0->0), fcn. (90->6), ass. (0->14)
t44 = qJ(2) + pkin(10);
t42 = sin(t44);
t43 = cos(t44);
t45 = sin(qJ(5));
t47 = cos(qJ(5));
t52 = -t42 * t47 + t43 * t45;
t49 = t42 * t45 + t43 * t47;
t48 = cos(qJ(1));
t46 = sin(qJ(1));
t38 = t49 * t48;
t37 = t52 * t48;
t36 = t49 * t46;
t35 = t52 * t46;
t1 = [-t36, t37, 0, 0, -t37, 0; t38, t35, 0, 0, -t35, 0; 0, t49, 0, 0, -t49, 0; t35, t38, 0, 0, -t38, 0; -t37, t36, 0, 0, -t36, 0; 0, -t52, 0, 0, t52, 0; -t48, 0, 0, 0, 0, 0; -t46, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
